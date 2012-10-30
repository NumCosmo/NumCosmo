/***************************************************************************
 *            ncm_model_ctrl.c
 *
 *  Mon February 27 12:10:09 2012
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
 * SECTION:ncm_model_ctrl
 * @title: Model update control object
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_model_ctrl.h"

enum
{
  PROP_0,
  PROP_MODEL,
};

G_DEFINE_TYPE (NcmModelCtrl, ncm_model_ctrl, G_TYPE_OBJECT);

static void
ncm_model_ctrl_init (NcmModelCtrl *ctrl)
{
  ctrl->model = NULL;
  ctrl->pkey = 0;
}

static void
ncm_model_ctrl_dispose (GObject *object)
{
  NcmModelCtrl *ctrl = NCM_MODEL_CTRL (object);
  if (ctrl->model)
	g_object_unref (ctrl->model);
  G_OBJECT_CLASS (ncm_model_ctrl_parent_class)->dispose (object);
}

static void
ncm_model_ctrl_finalize (GObject *object)
{
  G_OBJECT_CLASS (ncm_model_ctrl_parent_class)->finalize (object);
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
 * ncm_model_ctrl_copy:
 * @ctrl: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
   */
NcmModelCtrl *
ncm_model_ctrl_copy (NcmModelCtrl *ctrl)
{
  NcmModelCtrl *ctrl_copy = ncm_model_ctrl_new (ctrl->model);
  return ctrl_copy;
}

/**
 * ncm_model_ctrl_update:
 * @ctrl: FIXME
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */

/**
 * ncm_model_ctrl_get_model:
 * @ctrl: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
   */
NcmModel *
ncm_model_ctrl_get_model (NcmModelCtrl *ctrl)
{
  if (ctrl->model)
	g_object_ref (ctrl->model);
  return ctrl->model;
}

/**
 * ncm_model_ctrl_set_model:
 * @ctrl: FIXME
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_model_ctrl_set_model (NcmModelCtrl *ctrl, NcmModel *model)
{
  if (model != ctrl->model)
  {
	if (ctrl->model)
	  g_object_unref (ctrl->model);
	g_object_ref (model);
	ctrl->model = model;
	ctrl->pkey = model->pkey;
	return TRUE;
  }
  return FALSE;
}

/**
 * ncm_model_ctrl_has_model:
 * @ctrl: FIXME
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_model_ctrl_has_model (NcmModelCtrl *ctrl, NcmModel *model)
{
  return (model == ctrl->model);
}

/**
 * ncm_model_ctrl_has_model:
 * @ctrl: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
ncm_model_ctrl_force_update (NcmModelCtrl *ctrl)
{
  if (ctrl->model)
	g_object_unref (ctrl->model);
  ctrl->model = NULL;
  return;
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
	  g_value_set_object (value, ctrl->model);
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

  g_object_class_install_property (object_class,
                                   PROP_MODEL,
                                   g_param_spec_object ("model",
                                                        NULL,
                                                        "Last Model used",
                                                        NCM_TYPE_MODEL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  object_class->dispose = ncm_model_ctrl_dispose;
  object_class->finalize = ncm_model_ctrl_finalize;
}

void
ncm_model_ctrl_free (NcmModelCtrl *ctrl)
{
  g_object_unref (ctrl);
}

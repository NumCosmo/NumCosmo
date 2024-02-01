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
 * This object is employed to manage the status of a #NcmModel. It serves the purpose
 * of checking whether the model has been updated since the last call to
 * ncm_model_ctrl_update(). Calculation objects dependent on the model can utilize
 * this object to determine if updates are necessary.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_model_ctrl.h"
#include "math/ncm_cfg.h"

struct _NcmModelCtrl
{
  /*< private >*/
  GObject parent_instance;
  GWeakRef model_wr;
  gulong pkey;
  gboolean last_update;
  GPtrArray *submodel_ctrl;
  GArray *submodel_last_update;
};

enum
{
  PROP_0,
  PROP_MODEL,
};

G_DEFINE_TYPE (NcmModelCtrl, ncm_model_ctrl, G_TYPE_OBJECT)

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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_model_ctrl_class_init (NcmModelCtrlClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  /*GObjectClass* parent_class = G_OBJECT_CLASS (klass); */

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
 * @model: (allow-none): a #NcmModel or %NULL
 *
 * Creates a new #NcmModelCtrl object.
 *
 * Returns: (transfer full): a #NcmModelCtrl
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
 * ncm_model_ctrl_free:
 * @ctrl: a #NcmModelCtrl
 *
 * Decreases the reference count of @ctrl. If the reference count
 * reaches zero, @ctrl is freed.
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
 * Checks if *@ctrl is not %NULL, and if so, decreases the reference count
 * of @ctrl. If the reference count reaches zero, @ctrl is freed. The
 * pointer to @ctrl is set to %NULL.
 *
 */
void
ncm_model_ctrl_clear (NcmModelCtrl **ctrl)
{
  g_clear_object (ctrl);
}

/**
 * ncm_model_ctrl_update:
 * @ctrl: a #NcmModelCtrl
 * @model: a #NcmModel
 *
 * Compares the model inside @ctrl with @model and updates the status of
 * @ctrl. If the model inside @ctrl differs from @model, @ctrl is
 * updated, and TRUE is returned. Otherwise, FALSE is returned.
 *
 * If the model is the same but the model's pkey is different, @ctrl is
 * updated, and TRUE is returned. Otherwise, FALSE is returned.
 *
 * Submodels inside @model are also analyzed similarly, and @ctrl is
 * updated if necessary. If any submodel is updated, TRUE is returned.
 * Otherwise, FALSE is returned.
 *
 * Returns: TRUE if @ctrl was updated.
 */
gboolean
ncm_model_ctrl_update (NcmModelCtrl *ctrl, NcmModel *model)
{
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);
  guint64 pkey         = ncm_model_state_get_pkey (model);
  gboolean up          = FALSE;

  ctrl->last_update = FALSE;

  if (ctrl_model != model)
  {
    ncm_model_ctrl_set_model (ctrl, model);
    ctrl->last_update = TRUE;
  }
  else if (ctrl->pkey != pkey)
  {
    ctrl->pkey        = pkey;
    ctrl->last_update = TRUE;
  }

  up = up || ctrl->last_update;

  {
    const guint n = ncm_model_get_submodel_len (model);
    guint i;

    g_array_set_size (ctrl->submodel_last_update, n);

    for (i = 0; i < n; i++)
    {
      NcmModel *submodel = ncm_model_peek_submodel (model, i);

      g_array_index (ctrl->submodel_last_update, gboolean, i) = FALSE;

      if (i >= ctrl->submodel_ctrl->len)
      {
        NcmModelCtrl *sub_ctrl = ncm_model_ctrl_new (submodel);

        g_ptr_array_add (ctrl->submodel_ctrl, sub_ctrl);

        g_array_index (ctrl->submodel_last_update, gboolean, i) = TRUE;
      }
      else
      {
        NcmModelCtrl *sub_ctrl = g_ptr_array_index (ctrl->submodel_ctrl, i);

        if (ncm_model_ctrl_update (sub_ctrl, submodel))
          g_array_index (ctrl->submodel_last_update, gboolean, i) = TRUE;
      }

      up = up || g_array_index (ctrl->submodel_last_update, gboolean, i);
    }

    g_ptr_array_set_size (ctrl->submodel_ctrl, n);
  }

  ncm_model_clear (&ctrl_model);

  return up;
}

/**
 * ncm_model_ctrl_model_update:
 * @ctrl: a #NcmModelCtrl
 * @model: a #NcmModel
 *
 * Same as ncm_model_ctrl_update(), but only checks if the model objects
 * are the same. The pkey is not checked.
 *
 * Returns: TRUE if @ctrl was updated.
 */
gboolean
ncm_model_ctrl_model_update (NcmModelCtrl *ctrl, NcmModel *model)
{
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);
  gboolean up          = FALSE;

  ctrl->last_update = FALSE;

  if (ctrl_model != model)
  {
    ncm_model_ctrl_set_model (ctrl, model);
    ctrl->last_update = TRUE;
  }

  up = up || ctrl->last_update;

  {
    const guint n = ncm_model_get_submodel_len (model);
    guint i;

    g_array_set_size (ctrl->submodel_last_update, n);

    for (i = 0; i < n; i++)
    {
      NcmModel *submodel = ncm_model_peek_submodel (model, i);

      g_array_index (ctrl->submodel_last_update, gboolean, i) = FALSE;

      if (i >= ctrl->submodel_ctrl->len)
      {
        NcmModelCtrl *sub_ctrl = ncm_model_ctrl_new (submodel);

        g_ptr_array_add (ctrl->submodel_ctrl, sub_ctrl);

        g_array_index (ctrl->submodel_last_update, gboolean, i) = TRUE;
      }
      else
      {
        NcmModelCtrl *sub_ctrl = g_ptr_array_index (ctrl->submodel_ctrl, i);

        if (ncm_model_ctrl_model_update (sub_ctrl, submodel))
          g_array_index (ctrl->submodel_last_update, gboolean, i) = TRUE;
      }

      up = up || g_array_index (ctrl->submodel_last_update, gboolean, i);
    }

    g_ptr_array_set_size (ctrl->submodel_ctrl, n);
  }

  ncm_model_clear (&ctrl_model);

  return up;
}

/**
 * ncm_model_ctrl_get_model:
 * @ctrl: a #NcmModelCtrl
 *
 * Gets the current model inside @ctrl.
 *
 * Returns: (transfer full): a #NcmModel
 */
NcmModel *
ncm_model_ctrl_get_model (NcmModelCtrl *ctrl)
{
  return g_weak_ref_get (&ctrl->model_wr);
}

/**
 * ncm_model_ctrl_model_last_update:
 * @ctrl: a #NcmModelCtrl
 *
 * Checks if the main model was updated during the last call to
 * ncm_model_ctrl_update().
 *
 * Returns: TRUE if the main model was updated.
 */
gboolean
ncm_model_ctrl_model_last_update (NcmModelCtrl *ctrl)
{
  return ctrl->last_update;
}

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
gboolean
ncm_model_ctrl_model_has_submodel (NcmModelCtrl *ctrl, NcmModelID mid)
{
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);

  if (ctrl_model == NULL)
  {
    g_error ("ncm_model_ctrl_model_has_submodel: empty ctrl object.");

    return FALSE;
  }
  else
  {
    NcmModel *submodel    = ncm_model_peek_submodel_by_mid (ctrl_model, mid);
    gboolean has_submodel = FALSE;

    if (submodel != NULL)
      has_submodel = TRUE;

    ncm_model_clear (&ctrl_model);

    return has_submodel;
  }
}

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
gboolean
ncm_model_ctrl_submodel_last_update (NcmModelCtrl *ctrl, NcmModelID mid)
{
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);

  if (ctrl_model == NULL)
  {
    g_error ("ncm_model_ctrl_submodel_last_update: empty ctrl object.");

    return FALSE;
  }
  else
  {
    gint pos    = ncm_model_peek_submodel_pos_by_mid (ctrl_model, mid);
    gboolean up = FALSE;

    if (pos < 0)
      g_error ("ncm_model_ctrl_submodel_last_update: submodel `%s' not found in the main model `%s'.",
               ncm_mset_get_ns_by_id (mid),
               G_OBJECT_TYPE_NAME (ctrl_model));

    if ((guint) pos >= ctrl->submodel_last_update->len)
      g_error ("ncm_model_ctrl_submodel_last_update: submodel `%s' not found in ctrl object.\n"
               "ncm_model_ctrl_update() must always be called before ncm_model_ctrl_model_last_update_submodel().",
               ncm_mset_get_ns_by_id (mid));

    up = up || g_array_index (ctrl->submodel_last_update, gboolean, pos);

    ncm_model_clear (&ctrl_model);

    return up;
  }
}

/**
 * ncm_model_ctrl_set_model:
 * @ctrl: a #NcmModelCtrl
 * @model: a #NcmModel
 *
 * Sets the model inside @ctrl to @model.
 *
 * Returns: TRUE if @ctrl was updated.
 */
gboolean
ncm_model_ctrl_set_model (NcmModelCtrl *ctrl, NcmModel *model)
{
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);
  guint64 pkey         = ncm_model_state_get_pkey (model);
  gboolean up          = FALSE;

  if (model != ctrl_model)
  {
    g_weak_ref_set (&ctrl->model_wr, model);
    ctrl->pkey = pkey;
    up         = TRUE;
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
          up = TRUE;
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
 * @model: a #NcmModel
 *
 * Checks if the model inside @ctrl is the same as @model.
 *
 * Returns: TRUE if @model is the same as the model inside @ctrl.
 */
gboolean
ncm_model_ctrl_has_model (NcmModelCtrl *ctrl, NcmModel *model)
{
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);
  gboolean is_model    = (model == ctrl_model);

  ncm_model_clear (&ctrl_model);

  return is_model;
}

/**
 * ncm_model_ctrl_force_update:
 * @ctrl: a #NcmModelCtrl
 *
 * Forces an update on @ctrl. In practice, this function clears the
 * model inside @ctrl and all submodels.
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


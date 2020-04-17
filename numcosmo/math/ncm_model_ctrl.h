/***************************************************************************
 *            ncm_model_ctrl.h
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

#ifndef _NCM_MODEL_CTRL_H_
#define _NCM_MODEL_CTRL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NCM_TYPE_MODEL_CTRL             (ncm_model_ctrl_get_type ())
#define NCM_MODEL_CTRL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MODEL_CTRL, NcmModelCtrl))
#define NCM_MODEL_CTRL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MODEL_CTRL, NcmModelCtrlClass))
#define NCM_IS_MODEL_CTRL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MODEL_CTRL))
#define NCM_IS_MODEL_CTRL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MODEL_CTRL))
#define NCM_MODEL_CTRL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MODEL_CTRL, NcmModelCtrlClass))

typedef struct _NcmModelCtrlClass NcmModelCtrlClass;
typedef struct _NcmModelCtrl NcmModelCtrl;

struct _NcmModelCtrlClass
{
  /*< private >*/
  GObjectClass parent_class;
};

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

GType ncm_model_ctrl_get_type (void) G_GNUC_CONST;

NcmModelCtrl *ncm_model_ctrl_new (NcmModel *model);
gboolean ncm_model_ctrl_set_model (NcmModelCtrl *ctrl, NcmModel *model);
void ncm_model_ctrl_force_update (NcmModelCtrl *ctrl);
void ncm_model_ctrl_free (NcmModelCtrl *ctrl);
void ncm_model_ctrl_clear (NcmModelCtrl **ctrl);

NCM_INLINE NcmModel *ncm_model_ctrl_get_model (NcmModelCtrl *ctrl);
NCM_INLINE gboolean ncm_model_ctrl_update (NcmModelCtrl *ctrl, NcmModel *model);
NCM_INLINE gboolean ncm_model_ctrl_model_update (NcmModelCtrl *ctrl, NcmModel *model);

NCM_INLINE gboolean ncm_model_ctrl_model_last_update (NcmModelCtrl *ctrl);
NCM_INLINE gboolean ncm_model_ctrl_model_has_submodel (NcmModelCtrl *ctrl, NcmModelID mid);
NCM_INLINE gboolean ncm_model_ctrl_submodel_last_update (NcmModelCtrl *ctrl, NcmModelID mid);

G_END_DECLS

#endif /* _NCM_MODEL_CTRL_H_ */

#ifndef _NCM_MODEL_CTRL_INLINE_H_
#define _NCM_MODEL_CTRL_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE NcmModel *
ncm_model_ctrl_get_model (NcmModelCtrl *ctrl)
{
  return g_weak_ref_get (&ctrl->model_wr);
}

NCM_INLINE gboolean
ncm_model_ctrl_update (NcmModelCtrl *ctrl, NcmModel *model)
{
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);
  gboolean up = FALSE;

  ctrl->last_update = FALSE;
  if (ctrl_model != model)
  {
    ncm_model_ctrl_set_model (ctrl, model);
    ctrl->last_update = TRUE;
  }
  else if (ctrl->pkey != model->pkey)
  {
    ctrl->pkey = model->pkey;
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
        {
          g_array_index (ctrl->submodel_last_update, gboolean, i) = TRUE;
        }
      }
      up = up || g_array_index (ctrl->submodel_last_update, gboolean, i);
    }
    g_ptr_array_set_size (ctrl->submodel_ctrl, n);
  }

  ncm_model_clear (&ctrl_model);
  return up;
}

NCM_INLINE gboolean
ncm_model_ctrl_model_update (NcmModelCtrl *ctrl, NcmModel *model)
{
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);
  gboolean up = FALSE;

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
        {
          g_array_index (ctrl->submodel_last_update, gboolean, i) = TRUE;
        }
      }
      up = up || g_array_index (ctrl->submodel_last_update, gboolean, i);
    }
    g_ptr_array_set_size (ctrl->submodel_ctrl, n);
  }
  
  ncm_model_clear (&ctrl_model);  
  return up;
}

NCM_INLINE gboolean 
ncm_model_ctrl_model_last_update (NcmModelCtrl *ctrl)
{
  return ctrl->last_update;
}

NCM_INLINE gboolean 
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
    NcmModel *submodel = ncm_model_peek_submodel_by_mid (ctrl_model, mid);
    gboolean has_submodel = FALSE;

    if (submodel != NULL)
      has_submodel = TRUE;

    ncm_model_clear (&ctrl_model);

    return has_submodel;
  }
}

NCM_INLINE gboolean 
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
    gint pos = ncm_model_peek_submodel_pos_by_mid (ctrl_model, mid);
    gboolean up = FALSE;

    if (pos < 0)
      g_error ("ncm_model_ctrl_submodel_last_update: submodel `%s' not found in the main model `%s'.",
               ncm_mset_get_ns_by_id (mid),
               G_OBJECT_TYPE_NAME (ctrl_model));
    if (pos >= ctrl->submodel_last_update->len)
      g_error ("ncm_model_ctrl_submodel_last_update: submodel `%s' not found in ctrl object.\n" 
               "ncm_model_ctrl_update() must always be called before ncm_model_ctrl_model_last_update_submodel().",
               ncm_mset_get_ns_by_id (mid));

    up = up || g_array_index (ctrl->submodel_last_update, gboolean, pos);

    ncm_model_clear (&ctrl_model);

    return up;
  }
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_MODEL_CTRL_INLINE_H_ */

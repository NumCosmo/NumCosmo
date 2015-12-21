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

G_BEGIN_DECLS

#define NCM_TYPE_MODEL_CTRL             (ncm_model_ctrl_get_type ())
#define NCM_MODEL_CTRL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MODEL_CTRL, NcmModelCtrl))
#define NCM_MODEL_CTRL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MODEL_CTRL, NcmModelCtrlClass))
#define NCM_IS_MODEL_CTRL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MODEL_CTRL))
#define NCM_IS_MODEL_CTRL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MODEL_CTRL))
#define NCM_MODEL_CTRL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MODEL_CTRL, NcmModelCtrlClass))

typedef struct _NcmModelCtrlClass NcmModelCtrlClass;
typedef struct _NcmModelCtrl NcmModelCtrl;

struct _NcmModelCtrl
{
  /*< private >*/
  GObject parent_instance;
  NcmModel *model;
  gulong pkey;
};

struct _NcmModelCtrlClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_model_ctrl_get_type (void) G_GNUC_CONST;

NcmModelCtrl *ncm_model_ctrl_new (NcmModel *model);
NcmModelCtrl *ncm_model_ctrl_copy (NcmModelCtrl *ctrl);
gboolean ncm_model_ctrl_set_model (NcmModelCtrl *ctrl, NcmModel *model);
NcmModel *ncm_model_ctrl_get_model (NcmModelCtrl *ctrl);
void ncm_model_ctrl_force_update (NcmModelCtrl *ctrl);
void ncm_model_ctrl_free (NcmModelCtrl *ctrl);
void ncm_model_ctrl_clear (NcmModelCtrl **ctrl);

G_INLINE_FUNC NcmModel *ncm_model_ctrl_peek_model (NcmModelCtrl *ctrl);
G_INLINE_FUNC gboolean ncm_model_ctrl_update (NcmModelCtrl *ctrl, NcmModel *model);
G_INLINE_FUNC gboolean ncm_model_ctrl_model_update (NcmModelCtrl *ctrl, NcmModel *model);

G_END_DECLS

#endif /* _NCM_MODEL_CTRL_H_ */

#ifndef _NCM_MODEL_CTRL_INLINE_H_
#define _NCM_MODEL_CTRL_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC NcmModel *
ncm_model_ctrl_peek_model (NcmModelCtrl *ctrl)
{
  return ctrl->model;
}

G_INLINE_FUNC gboolean
ncm_model_ctrl_update (NcmModelCtrl *ctrl, NcmModel *model)
{
  if (ctrl->model != model)
  {
    ncm_model_ctrl_set_model (ctrl, model);
    return TRUE;
  }
  else if (ctrl->pkey != model->pkey)
  {
    ctrl->pkey = model->pkey;
    return TRUE;
  }
  return FALSE;
}

G_INLINE_FUNC gboolean
ncm_model_ctrl_model_update (NcmModelCtrl *ctrl, NcmModel *model)
{
  if (ctrl->model != model)
  {
    ncm_model_ctrl_set_model (ctrl, model);
    return TRUE;
  }
  return FALSE;
}


G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_MODEL_CTRL_INLINE_H_ */

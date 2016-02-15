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

#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32))
#define NCM_MODEL_CTRL_USE_WEAKREF 1
#endif

struct _NcmModelCtrl
{
  /*< private >*/
  GObject parent_instance;
#ifdef NCM_MODEL_CTRL_USE_WEAKREF
  GWeakRef model_wf;
#else /* NCM_MODEL_CTRL_USE_WEAKREF */
  NcmModel *model;
#endif  /* NCM_MODEL_CTRL_USE_WEAKREF */
  gulong pkey;
};

struct _NcmModelCtrlClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_model_ctrl_get_type (void) G_GNUC_CONST;

NcmModelCtrl *ncm_model_ctrl_new (NcmModel *model);
gboolean ncm_model_ctrl_set_model (NcmModelCtrl *ctrl, NcmModel *model);
void ncm_model_ctrl_force_update (NcmModelCtrl *ctrl);
void ncm_model_ctrl_free (NcmModelCtrl *ctrl);
void ncm_model_ctrl_clear (NcmModelCtrl **ctrl);

G_INLINE_FUNC NcmModel *ncm_model_ctrl_get_model (NcmModelCtrl *ctrl);
G_INLINE_FUNC gboolean ncm_model_ctrl_update (NcmModelCtrl *ctrl, NcmModel *model);
G_INLINE_FUNC gboolean ncm_model_ctrl_model_update (NcmModelCtrl *ctrl, NcmModel *model);

G_END_DECLS

#endif /* _NCM_MODEL_CTRL_H_ */

#ifndef _NCM_MODEL_CTRL_INLINE_H_
#define _NCM_MODEL_CTRL_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC NcmModel *
ncm_model_ctrl_get_model (NcmModelCtrl *ctrl)
{
#ifdef NCM_MODEL_CTRL_USE_WEAKREF
  return g_weak_ref_get (&ctrl->model_wf);
#else
  return ncm_model_ref (ctrl->model);
#endif
}

G_INLINE_FUNC gboolean
ncm_model_ctrl_update (NcmModelCtrl *ctrl, NcmModel *model)
{
  gboolean up = FALSE;
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);
  if (ctrl_model != model)
  {
    ncm_model_ctrl_set_model (ctrl, model);
    up = TRUE;
  }
  else if (ctrl->pkey != model->pkey)
  {
    ctrl->pkey = model->pkey;
    up = TRUE;
  }
  ncm_model_clear (&ctrl_model);  
  return up;
}

G_INLINE_FUNC gboolean
ncm_model_ctrl_model_update (NcmModelCtrl *ctrl, NcmModel *model)
{
  gboolean up = FALSE;
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);
  if (ctrl_model != model)
  {
    ncm_model_ctrl_set_model (ctrl, model);
    up = TRUE;
  }
  ncm_model_clear (&ctrl_model);  
  return up;
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_MODEL_CTRL_INLINE_H_ */

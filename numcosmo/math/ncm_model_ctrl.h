/***************************************************************************
 *            ncm_model_ctrl.h
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

#ifndef _NCM_MODEL_CTRL_H_
#define _NCM_MODEL_CTRL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NCM_TYPE_MODEL_CTRL (ncm_model_ctrl_get_type ())

G_DECLARE_FINAL_TYPE (NcmModelCtrl, ncm_model_ctrl, NCM, MODEL_CTRL, GObject)

NcmModelCtrl *ncm_model_ctrl_new (NcmModel * model);
gboolean ncm_model_ctrl_set_model (NcmModelCtrl *ctrl, NcmModel *model);
void ncm_model_ctrl_force_update (NcmModelCtrl *ctrl);
void ncm_model_ctrl_free (NcmModelCtrl *ctrl);
void ncm_model_ctrl_clear (NcmModelCtrl **ctrl);

NcmModel *ncm_model_ctrl_get_model (NcmModelCtrl *ctrl);
gboolean ncm_model_ctrl_update (NcmModelCtrl *ctrl, NcmModel *model);
gboolean ncm_model_ctrl_model_update (NcmModelCtrl *ctrl, NcmModel *model);

gboolean ncm_model_ctrl_model_last_update (NcmModelCtrl *ctrl);
gboolean ncm_model_ctrl_model_has_submodel (NcmModelCtrl *ctrl, NcmModelID mid);
gboolean ncm_model_ctrl_submodel_last_update (NcmModelCtrl *ctrl, NcmModelID mid);

G_END_DECLS

#endif /* _NCM_MODEL_CTRL_H_ */


/***************************************************************************
 *            nc_hireion_camb_reparam_tau.h
 *
 *  Tue December 15 04:31:46 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hireion_camb_reparam_tau.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_HIREION_CAMB_REPARAM_TAU_H_
#define _NC_HIREION_CAMB_REPARAM_TAU_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_reparam.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HIREION_CAMB_REPARAM_TAU             (nc_hireion_camb_reparam_tau_get_type ())
#define NC_HIREION_CAMB_REPARAM_TAU(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIREION_CAMB_REPARAM_TAU, NcHIReionCambReparamTau))
#define NC_HIREION_CAMB_REPARAM_TAU_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIREION_CAMB_REPARAM_TAU, NcHIReionCambReparamTauClass))
#define NC_IS_HIREION_CAMB_REPARAM_TAU(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIREION_CAMB_REPARAM_TAU))
#define NC_IS_HIREION_CAMB_REPARAM_TAU_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIREION_CAMB_REPARAM_TAU))
#define NC_HIREION_CAMB_REPARAM_TAU_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIREION_CAMB_REPARAM_TAU, NcHIReionCambReparamTauClass))

typedef struct _NcHIReionCambReparamTauClass NcHIReionCambReparamTauClass;
typedef struct _NcHIReionCambReparamTau NcHIReionCambReparamTau;

struct _NcHIReionCambReparamTauClass
{
  /*< private >*/
  NcmReparamClass parent_class;
};

struct _NcHIReionCambReparamTau
{
  /*< private >*/
  NcmReparam parent_instance;
  NcmModelCtrl *ctrl;
};

GType nc_hireion_camb_reparam_tau_get_type (void) G_GNUC_CONST;

NcHIReionCambReparamTau *nc_hireion_camb_reparam_tau_new (guint length, NcHICosmo *cosmo);

G_END_DECLS

#endif /* _NC_HIREION_CAMB_REPARAM_TAU_H_ */

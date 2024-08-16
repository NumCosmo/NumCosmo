/***************************************************************************
 *            nc_hicosmo_de_reparam_cmb.h
 *
 *  Fri April 15 16:14:17 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hicosmo_de_reparam_cmb.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_HICOSMO_DE_REPARAM_CMB_H_
#define _NC_HICOSMO_DE_REPARAM_CMB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_reparam.h>
#include <numcosmo/math/ncm_model.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE_REPARAM_CMB             (nc_hicosmo_de_reparam_cmb_get_type ())
#define NC_HICOSMO_DE_REPARAM_CMB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_DE_REPARAM_CMB, NcHICosmoDEReparamCMB))
#define NC_HICOSMO_DE_REPARAM_CMB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_DE_REPARAM_CMB, NcHICosmoDEReparamCMBClass))
#define NC_IS_HICOSMO_DE_REPARAM_CMB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_DE_REPARAM_CMB))
#define NC_IS_HICOSMO_DE_REPARAM_CMB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_DE_REPARAM_CMB))
#define NC_HICOSMO_DE_REPARAM_CMB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_DE_REPARAM_CMB, NcHICosmoDEReparamCMBClass))

typedef struct _NcHICosmoDEReparamCMBClass NcHICosmoDEReparamCMBClass;
typedef struct _NcHICosmoDEReparamCMB NcHICosmoDEReparamCMB;

struct _NcHICosmoDEReparamCMBClass
{
  /*< private >*/
  NcmReparamClass parent_class;
};

struct _NcHICosmoDEReparamCMB
{
  /*< private >*/
  NcmReparam parent_instance;
};

GType nc_hicosmo_de_reparam_cmb_get_type (void) G_GNUC_CONST;

NcHICosmoDEReparamCMB *nc_hicosmo_de_reparam_cmb_new (guint length);

G_END_DECLS

#endif /* _NC_HICOSMO_DE_REPARAM_CMB_H_ */


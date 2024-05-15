/***************************************************************************
 *            nc_hicosmo_de_reparam_ok.h
 *
 *  Mon October 26 10:50:22 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hicosmo_de_reparam_ok.h
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

#ifndef _NC_HICOSMO_DE_REPARAM_OK_H_
#define _NC_HICOSMO_DE_REPARAM_OK_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_reparam.h>
#include <numcosmo/math/ncm_model.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE_REPARAM_OK             (nc_hicosmo_de_reparam_ok_get_type ())
#define NC_HICOSMO_DE_REPARAM_OK(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_DE_REPARAM_OK, NcHICosmoDEReparamOk))
#define NC_HICOSMO_DE_REPARAM_OK_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_DE_REPARAM_OK, NcHICosmoDEReparamOkClass))
#define NC_IS_HICOSMO_DE_REPARAM_OK(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_DE_REPARAM_OK))
#define NC_IS_HICOSMO_DE_REPARAM_OK_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_DE_REPARAM_OK))
#define NC_HICOSMO_DE_REPARAM_OK_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_DE_REPARAM_OK, NcHICosmoDEReparamOkClass))

typedef struct _NcHICosmoDEReparamOkClass NcHICosmoDEReparamOkClass;
typedef struct _NcHICosmoDEReparamOk NcHICosmoDEReparamOk;

struct _NcHICosmoDEReparamOkClass
{
  /*< private >*/
  NcmReparamClass parent_class;
};

struct _NcHICosmoDEReparamOk
{
  /*< private >*/
  NcmReparam parent_instance;
};

GType nc_hicosmo_de_reparam_ok_get_type (void) G_GNUC_CONST;

NcHICosmoDEReparamOk *nc_hicosmo_de_reparam_ok_new (guint length);

G_END_DECLS

#endif /* _NC_HICOSMO_DE_REPARAM_OK_H_ */


/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_de_cont.h
 *
 *  Thu December 15 15:08:26 2020
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_de_cont.h
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_DE_CONT_H_
#define _NC_DE_CONT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_csq1d.h>
#include <numcosmo/math/ncm_spline.h>

G_BEGIN_DECLS

#define NC_TYPE_DE_CONT             (nc_de_cont_get_type ())
#define NC_DE_CONT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DE_CONT, NcDECont))
#define NC_DE_CONT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DE_CONT, NcDEContClass))
#define NC_IS_DE_CONT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DE_CONT))
#define NC_IS_DE_CONT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DE_CONT))
#define NC_DE_CONT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DE_CONT, NcDEContClass))

typedef struct _NcDEContClass NcDEContClass;
typedef struct _NcDECont NcDECont;
typedef struct _NcDEContPrivate NcDEContPrivate;

struct _NcDEContClass
{
  /*< private >*/
  NcmCSQ1DClass parent_class;
};

struct _NcDECont
{
  /*< private >*/
  NcmCSQ1D parent_instance;
  NcDEContPrivate *priv;
};

GType nc_de_cont_get_type (void) G_GNUC_CONST;

NcDECont *nc_de_cont_new (const gdouble Omegaw, const gdouble OmegaL, const gdouble cs2, const gdouble w);
NcDECont *nc_de_cont_ref (NcDECont *dec);

void nc_de_cont_free (NcDECont *dec);
void nc_de_cont_clear (NcDECont **dec);

G_END_DECLS

#endif /* _NC_DE_CONT_H_ */

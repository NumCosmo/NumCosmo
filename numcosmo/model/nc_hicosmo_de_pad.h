/***************************************************************************
 *            nc_hicosmo_de_pad.h
 *
 *  Mon Aug 11 19:57:55 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NC_HICOSMO_DE_PAD_H_
#define _NC_HICOSMO_DE_PAD_H_

#include <glib-object.h>
#include <numcosmo/model/nc_hicosmo_de.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE_PAD             (nc_hicosmo_de_pad_get_type ())
#define NC_HICOSMO_DE_PAD(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_DE_PAD, NcHICosmoDEPad))
#define NC_HICOSMO_DE_PAD_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_DE_PAD, NcHICosmoDEPadClass))
#define NC_IS_HICOSMO_DE_PAD(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_DE_PAD))
#define NC_IS_HICOSMO_DE_PAD_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_DE_PAD))
#define NC_HICOSMO_DE_PAD_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_DE_PAD, NcHICosmoDEPadClass))

typedef struct _NcHICosmoDEPadClass NcHICosmoDEPadClass;
typedef struct _NcHICosmoDEPad NcHICosmoDEPad;


/**
 * NcHICosmoDEPadParams:
 * @NC_HICOSMO_DE_PAD_W0: FIXME
 * @NC_HICOSMO_DE_PAD_W1: FIXME
 *
 * FIXME
 */
typedef enum _NcHICosmoDEPadParams
{
  NC_HICOSMO_DE_PAD_W0 = NC_HICOSMO_DE_SPARAM_LEN,
  NC_HICOSMO_DE_PAD_W1,         /*< private >*/
  NC_HICOSMO_DE_PAD_SPARAM_LEN, /*< skip >*/
} NcHICosmoDEPadParams;

#define NC_HICOSMO_DE_PAD_DEFAULT_W0 (-1.0)
#define NC_HICOSMO_DE_PAD_DEFAULT_W1 ( 0.0)

#define NC_HICOSMO_DE_PM_N (NC_HICOSMO_DE_PAD_W1 + 1 - NC_HICOSMO_DE_BASE_N

struct _NcHICosmoDEPadClass
{
	/*< private >*/
	NcHICosmoDEClass parent_class;
};

struct _NcHICosmoDEPad
{
	/*< private >*/
	NcHICosmoDE parent_instance;
};

GType nc_hicosmo_de_pad_get_type (void) G_GNUC_CONST;

NcHICosmoDEPad *nc_hicosmo_de_pad_new (void);

G_END_DECLS

#endif /* _NC_HICOSMO_DE_PAD_H_ */

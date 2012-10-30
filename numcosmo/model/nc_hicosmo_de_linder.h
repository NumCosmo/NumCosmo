/***************************************************************************
 *            nc_hicosmo_de_linder.h
 *
 *  Mon Aug 11 19:56:54 2008
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

#ifndef _NC_HICOSMO_DE_LINDER_H_
#define _NC_HICOSMO_DE_LINDER_H_

#include <glib-object.h>
#include <numcosmo/model/nc_hicosmo_de.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE_LINDER             (nc_hicosmo_de_linder_get_type ())
#define NC_HICOSMO_DE_LINDER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_DE_LINDER, NcHICosmoDELinder))
#define NC_HICOSMO_DE_LINDER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_DE_LINDER, NcHICosmoDELinderClass))
#define NC_IS_HICOSMO_DE_LINDER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_DE_LINDER))
#define NC_IS_HICOSMO_DE_LINDER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_DE_LINDER))
#define NC_HICOSMO_DE_LINDER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_DE_LINDER, NcHICosmoDELinderClass))

typedef struct _NcHICosmoDELinderClass NcHICosmoDELinderClass;
typedef struct _NcHICosmoDELinder NcHICosmoDELinder;

/**
 * NcHICosmoDELinderParams:
 * @NC_HICOSMO_DE_LINDER_W0: FIXME
 * @NC_HICOSMO_DE_LINDER_W1: FIXME
 *
 * FIXME
 */
typedef enum _NcHICosmoDELinderParams
{
  NC_HICOSMO_DE_LINDER_W0 = NC_HICOSMO_DE_SPARAM_LEN,
  NC_HICOSMO_DE_LINDER_W1,         /*< private >*/
  NC_HICOSMO_DE_LINDER_SPARAM_LEN, /*< skip >*/
} NcHICosmoDELinderParams;

#define NC_HICOSMO_DE_LINDER_DEFAULT_W0 (-1.0)
#define NC_HICOSMO_DE_LINDER_DEFAULT_W1 ( 0.0)

#define NC_HICOSMO_DE_LINDER_N (NC_HICOSMO_DE_LINDER_W1 + 1 - NC_HICOSMO_DE_BASE_N)

struct _NcHICosmoDELinderClass
{
	/*< private >*/
	NcHICosmoDEClass parent_class;
};

struct _NcHICosmoDELinder
{
	/*< private >*/
	NcHICosmoDE parent_instance;
};

GType nc_hicosmo_de_linder_get_type (void) G_GNUC_CONST;

NcHICosmoDELinder *nc_hicosmo_de_linder_new (void);

G_END_DECLS

#endif /* _NC_HICOSMO_DE_LINDER_H_ */

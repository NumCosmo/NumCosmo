/***************************************************************************
 *            nc_hicosmo_de_qe.h
 *
 *  Mon Aug 11 19:58:55 2008
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

#ifndef _NC_HICOSMO_DE_QE_H_
#define _NC_HICOSMO_DE_QE_H_

#include <glib-object.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE_QE             (nc_hicosmo_de_qe_get_type ())
#define NC_HICOSMO_DE_QE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_DE_QE, NcHICosmoDEQe))
#define NC_HICOSMO_DE_QE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_DE_QE, NcHICosmoDEQeClass))
#define NC_IS_HICOSMO_DE_QE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_DE_QE))
#define NC_IS_HICOSMO_DE_QE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_DE_QE))
#define NC_HICOSMO_DE_QE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_DE_QE, NcHICosmoDEQeClass))

typedef struct _NcHICosmoDEQeClass NcHICosmoDEQeClass;
typedef struct _NcHICosmoDEQe NcHICosmoDEQe;

/**
 * NcHICosmoDEQEParams:
 * @NC_HICOSMO_DE_QE_W0: FIXME
 * @NC_HICOSMO_DE_QE_W1: FIXME
 *
 * FIXME
 */
typedef enum _NcHICosmoDEQEParams
{
  NC_HICOSMO_DE_QE_W0 = NC_HICOSMO_DE_SPARAM_LEN,
  NC_HICOSMO_DE_QE_W1,         /*< private >*/
  NC_HICOSMO_DE_QE_SPARAM_LEN, /*< skip >*/
} NcHICosmoDEQEParams;

#define NC_HICOSMO_DE_QE_DEFAULT_W0 (-1.0)
#define NC_HICOSMO_DE_QE_DEFAULT_W1 ( 0.0)

#define NC_HICOSMO_DE_QE_N (NC_HICOSMO_DE_QE_W1 + 1 - NC_HICOSMO_DE_BASE_N)

struct _NcHICosmoDEQeClass
{
	/*< private >*/
	NcHICosmoDEClass parent_class;
};

struct _NcHICosmoDEQe
{
	/*< private >*/
	NcHICosmoDE parent_instance;
};

GType nc_hicosmo_de_qe_get_type (void) G_GNUC_CONST;

NcHICosmoDEQe *nc_hicosmo_de_qe_new (void);

G_END_DECLS

#endif /* _NC_HICOSMO_DE_QE_H_ */

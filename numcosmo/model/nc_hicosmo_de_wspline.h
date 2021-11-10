/***************************************************************************
 *            nc_hicosmo_de_wspline.c
 *
 *  Mon Oct 11 16:22:12 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 *  Copyright  2021  Sanderson Carlos Ribeiro
 *  <sander23.ribeiro@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2021 <sandro@isoftware.com.br>
 * Copyright (C) Sanderson Carlos Ribeiro 2021 <sander23.ribeiro@uel.br>
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

#ifndef _NC_HICOSMO_DE_WSPLINE_H_
#define _NC_HICOSMO_DE_WSPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/model/nc_hicosmo_de.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE_WSPLINE             (nc_hicosmo_de_wspline_get_type ())
#define NC_HICOSMO_DE_WSPLINE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_DE_WSPLINE, NcHICosmoDEWSpline))
#define NC_HICOSMO_DE_WSPLINE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_DE_WSPLINE, NcHICosmoDEWSplineClass))
#define NC_IS_HICOSMO_DE_WSPLINE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_DE_WSPLINE))
#define NC_IS_HICOSMO_DE_WSPLINE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_DE_WSPLINE))
#define NC_HICOSMO_DE_WSPLINE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_DE_WSPLINE, NcHICosmoDEWSplineClass))

typedef struct _NcHICosmoDEWSplineClass NcHICosmoDEWSplineClass;
typedef struct _NcHICosmoDEWSpline NcHICosmoDEWSpline;
typedef struct _NcHICosmoDEWSplinePrivate NcHICosmoDEWSplinePrivate;

/**
 * NcHICosmoDEWSplineSParams:
 * @NC_HICOSMO_DE_WSPLINE_W: constant parameter
 *
 * Dark Energy equation of state: $w(z) = w$.
 * 
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_DE_WSPLINE_SPARAMS >*/
{
  /*NC_HICOSMO_DE_WSPLINE_W = NC_HICOSMO_DE_SPARAM_LEN,*/
  /* < private > */
  NC_HICOSMO_DE_WSPLINE_SPARAM_LEN = NC_HICOSMO_DE_SPARAM_LEN,                   /*< skip >*/
} NcHICosmoDEWSplineSParams;


/**
 * NcHICosmoDEWSplineVParams:
 * @NC_HICOSMO_DE_WSPLINE_W: FIXME
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_DE_WSPLINE_VPARAMS >*/
{
  NC_HICOSMO_DE_WSPLINE_W = NC_HICOSMO_DE_VPARAM_LEN,
  /* < private > */
  NC_HICOSMO_DE_WSPLINE_VPARAM_LEN, /*< skip >*/
} NcHICosmoDEWSplineVParams;

#define NC_HICOSMO_DE_WSPLINE_DEFAULT_W0 (-1.0)

#define NC_HICOSMO_DE_WSPLINE_N (NC_HICOSMO_DE_WSPLINE_W + 1 - NC_HICOSMO_DE_BASE_N)

struct _NcHICosmoDEWSplineClass
{
  /*< private >*/
  NcHICosmoDEClass parent_class;
};

struct _NcHICosmoDEWSpline
{
  /*< private >*/
  NcHICosmoDE parent_instance;
  NcHICosmoDEWSplinePrivate *priv;
};

GType nc_hicosmo_de_wspline_get_type (void) G_GNUC_CONST;

NcHICosmoDEWSpline *nc_hicosmo_de_wspline_new (gsize nknots, const gdouble z_f);

NcmVector *nc_hicosmo_de_wspline_get_alpha (NcHICosmoDEWSpline *wspline);

G_END_DECLS

#endif /* _NC_HICOSMO_DE_WSPLINE_H_ */

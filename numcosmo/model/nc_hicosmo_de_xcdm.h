/***************************************************************************
 *            nc_hicosmo_de_xcdm.h
 *
 *  Mon Aug 11 19:57:16 2008
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

#ifndef _NC_HICOSMO_DE_XCDM_H_
#define _NC_HICOSMO_DE_XCDM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/model/nc_hicosmo_de.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE_XCDM             (nc_hicosmo_de_xcdm_get_type ())
#define NC_HICOSMO_DE_XCDM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_DE_XCDM, NcHICosmoDEXcdm))
#define NC_HICOSMO_DE_XCDM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_DE_XCDM, NcHICosmoDEXcdmClass))
#define NC_IS_HICOSMO_DE_XCDM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_DE_XCDM))
#define NC_IS_HICOSMO_DE_XCDM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_DE_XCDM))
#define NC_HICOSMO_DE_XCDM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_DE_XCDM, NcHICosmoDEXcdmClass))

typedef struct _NcHICosmoDEXcdmClass NcHICosmoDEXcdmClass;
typedef struct _NcHICosmoDEXcdm NcHICosmoDEXcdm;

/**
 * NcHICosmoDEXCDMSParams:
 * @NC_HICOSMO_DE_XCDM_W: constant parameter
 *
 * Dark Energy equation of state: $w(z) = w$.
 * 
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_DE_XCDM_SPARAMS >*/
{
  NC_HICOSMO_DE_XCDM_W = NC_HICOSMO_DE_SPARAM_LEN, 
  /* < private > */
  NC_HICOSMO_DE_XCDM_SPARAM_LEN,                   /*< skip >*/
} NcHICosmoDEXCDMSParams;

#define NC_HICOSMO_DE_XCDM_DEFAULT_W0 (-1.0)

#define NC_HICOSMO_DE_XCDM_N (NC_HICOSMO_DE_XCDM_W + 1 - NC_HICOSMO_DE_BASE_N)

struct _NcHICosmoDEXcdmClass
{
  /*< private >*/
  NcHICosmoDEClass parent_class;
};

struct _NcHICosmoDEXcdm
{
  /*< private >*/
  NcHICosmoDE parent_instance;
};

GType nc_hicosmo_de_xcdm_get_type (void) G_GNUC_CONST;

NcHICosmoDEXcdm *nc_hicosmo_de_xcdm_new (void);

G_END_DECLS

#endif /* _NC_HICOSMO_DE_XCDM_H_ */

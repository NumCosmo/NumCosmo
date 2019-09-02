/***************************************************************************
 *            nc_hicosmo_de_cpl.h
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

#ifndef _NC_HICOSMO_DE_CPL_H_
#define _NC_HICOSMO_DE_CPL_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/model/nc_hicosmo_de.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE_CPL             (nc_hicosmo_de_cpl_get_type ())
#define NC_HICOSMO_DE_CPL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_DE_CPL, NcHICosmoDECpl))
#define NC_HICOSMO_DE_CPL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_DE_CPL, NcHICosmoDECplClass))
#define NC_IS_HICOSMO_DE_CPL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_DE_CPL))
#define NC_IS_HICOSMO_DE_CPL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_DE_CPL))
#define NC_HICOSMO_DE_CPL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_DE_CPL, NcHICosmoDECplClass))

typedef struct _NcHICosmoDECplClass NcHICosmoDECplClass;
typedef struct _NcHICosmoDECpl NcHICosmoDECpl;

/**
 * NcHICosmoDECplSParams:
 * @NC_HICOSMO_DE_CPL_W0: constant parameter
 * @NC_HICOSMO_DE_CPL_W1: constant parameter which multiplies the redshift term
 *
 * Dark Energy equation of state: $w(z) = w_0 + w_1 \frac{z}{1.0 + z}$
 * 
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_DE_CPL_SPARAMS >*/
{
  NC_HICOSMO_DE_CPL_W0 = NC_HICOSMO_DE_SPARAM_LEN,
  NC_HICOSMO_DE_CPL_W1,         
  /* < private > */
  NC_HICOSMO_DE_CPL_SPARAM_LEN, /*< skip >*/
} NcHICosmoDECplSParams;

#define NC_HICOSMO_DE_CPL_DEFAULT_W0 (-1.0)
#define NC_HICOSMO_DE_CPL_DEFAULT_W1 ( 0.0)

#define NC_HICOSMO_DE_CPL_N (NC_HICOSMO_DE_CPL_W1 + 1 - NC_HICOSMO_DE_BASE_N)

struct _NcHICosmoDECplClass
{
  /*< private >*/
  NcHICosmoDEClass parent_class;
};

struct _NcHICosmoDECpl
{
  /*< private >*/
  NcHICosmoDE parent_instance;
};

GType nc_hicosmo_de_cpl_get_type (void) G_GNUC_CONST;

NcHICosmoDECpl *nc_hicosmo_de_cpl_new (void);

#ifdef HAVE_CCL
typedef struct ccl_parameters ccl_parameters;
NcHICosmoDECpl *nc_hicosmo_de_cpl_new_from_ccl (ccl_parameters *ccl_params);
#endif /* HAVE_CCL */

G_END_DECLS

#endif /* _NC_HICOSMO_DE_CPL_H_ */

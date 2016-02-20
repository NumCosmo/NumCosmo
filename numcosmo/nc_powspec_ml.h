/***************************************************************************
 *            nc_powspec_ml.h
 *
 *  Thu February 18 12:32:25 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_powspec_ml.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_POWSPEC_ML_H_
#define _NC_POWSPEC_ML_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_powspec.h>

G_BEGIN_DECLS

#define NC_TYPE_POWSPEC_ML             (nc_powspec_ml_get_type ())
#define NC_POWSPEC_ML(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_POWSPEC_ML, NcPowspecML))
#define NC_POWSPEC_ML_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_POWSPEC_ML, NcPowspecMLClass))
#define NC_IS_POWSPEC_ML(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_POWSPEC_ML))
#define NC_IS_POWSPEC_ML_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_POWSPEC_ML))
#define NC_POWSPEC_ML_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_POWSPEC_ML, NcPowspecMLClass))

typedef struct _NcPowspecMLClass NcPowspecMLClass;
typedef struct _NcPowspecML NcPowspecML;

struct _NcPowspecMLClass
{
  /*< private > */
  NcmPowspecClass parent_class;
};

struct _NcPowspecML
{
  /*< private > */
  NcmPowspec parent_instance;
};

GType nc_powspec_ml_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_POWSPEC_ML_H_ */

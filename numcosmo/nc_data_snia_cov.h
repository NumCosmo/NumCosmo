/***************************************************************************
 *            nc_data_snia_cov.h
 *
 *  Sat December 08 15:58:20 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_SNIA_COV_H_
#define _NC_DATA_SNIA_COV_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>
#include <numcosmo/nc_snia_dist_cov.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_SNIA_COV             (nc_data_snia_cov_get_type ())
#define NC_DATA_SNIA_COV(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_SNIA_COV, NcDataSNIACov))
#define NC_DATA_SNIA_COV_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_SNIA_COV, NcDataSNIACovClass))
#define NC_IS_DATA_SNIA_COV(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_SNIA_COV))
#define NC_IS_DATA_SNIA_COV_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_SNIA_COV))
#define NC_DATA_SNIA_COV_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_SNIA_COV, NcDataSNIACovClass))

typedef struct _NcDataSNIACovClass NcDataSNIACovClass;
typedef struct _NcDataSNIACov NcDataSNIACov;

struct _NcDataSNIACovClass
{
  /*< private >*/
  NcmDataGaussCovClass parent_class;
};

struct _NcDataSNIACov
{
  /*< private >*/
  NcmDataGaussCov parent_instance;
  NcSNIADistCov *dcov;
};

GType nc_data_snia_cov_get_type (void) G_GNUC_CONST;

NcmData *nc_data_snia_cov_new (NcSNIADistCov *dcov, gboolean use_det);

void nc_data_snia_cov_set_dcov (NcDataSNIACov *snia_cov, NcSNIADistCov *dcov);

G_END_DECLS

#endif /* _NC_DATA_SNIA_COV_H_ */


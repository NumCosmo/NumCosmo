/***************************************************************************
 *            nc_data_cmb_dist_priors.h
 *
 *  Thu Apr 22 15:56:59 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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

#ifndef _NC_DATA_CMB_DIST_PRIORS_H_
#define _NC_DATA_CMB_DIST_PRIORS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/nc_data_cmb.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CMB_DIST_PRIORS             (nc_data_cmb_dist_priors_get_type ())
#define NC_DATA_CMB_DIST_PRIORS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_CMB_DIST_PRIORS, NcDataCMBDistPriors))
#define NC_DATA_CMB_DIST_PRIORS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_CMB_DIST_PRIORS, NcDataCMBDistPriorsClass))
#define NC_IS_DATA_CMB_DIST_PRIORS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_CMB_DIST_PRIORS))
#define NC_IS_DATA_CMB_DIST_PRIORS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_CMB_DIST_PRIORS))
#define NC_DATA_CMB_DIST_PRIORS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_CMB_DIST_PRIORS, NcDataCMBDistPriorsClass))

typedef struct _NcDataCMBDistPriorsClass NcDataCMBDistPriorsClass;
typedef struct _NcDataCMBDistPriors NcDataCMBDistPriors;

struct _NcDataCMBDistPriors
{
  /*< private >*/
  NcmDataGauss parent_instance;
  NcDistance *dist;
  NcDataCMBId id;
};

struct _NcDataCMBDistPriorsClass
{
  /*< private >*/
  NcmDataGaussClass parent_class;
};

GType nc_data_cmb_dist_priors_get_type (void) G_GNUC_CONST;

NcmData *nc_data_cmb_dist_priors_new (NcDistance *dist, NcDataCMBId id);

void nc_data_cmb_dist_priors_set_sample (NcDataCMBDistPriors *cmb_dist_priors, NcDataCMBId id);
NcDataCMBId nc_data_cmb_dist_priors_get_sample (NcDataCMBDistPriors *cmb_dist_priors);

G_END_DECLS

#endif /* _NC_DATA_CMB_DIST_PRIORS_H_ */


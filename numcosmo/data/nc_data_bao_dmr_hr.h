/***************************************************************************
 *            nc_data_bao_dmr_hr.h
 *
 *  Mon Oct 10 23:29:40 2016
 *  Copyright  2016  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_bao_dmr_hr.h
 * Copyright (C) 2016 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_DATA_BAO_DMR_HR_H_
#define _NC_DATA_BAO_DMR_HR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/data/nc_data_bao.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_BAO_DMR_HR             (nc_data_bao_dmr_hr_get_type ())
#define NC_DATA_BAO_DMR_HR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_BAO_DMR_HR, NcDataBaoDMrHr))
#define NC_DATA_BAO_DMR_HR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_BAO_DMR_HR, NcDataBaoDMrHrClass))
#define NC_IS_DATA_BAO_DMR_HR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_BAO_DMR_HR))
#define NC_IS_DATA_BAO_DMR_HR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_BAO_DMR_HR))
#define NC_DATA_BAO_DMR_HR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_BAO_DMR_HR, NcDataBaoDMrHrClass))

typedef struct _NcDataBaoDMrHrClass NcDataBaoDMrHrClass;
typedef struct _NcDataBaoDMrHr NcDataBaoDMrHr;

struct _NcDataBaoDMrHrClass
{
  /*< private >*/
  NcmDataGaussCovClass parent_class;
};

struct _NcDataBaoDMrHr
{
  /*< private >*/
  NcmDataGaussCov parent_instance;
  NcDistance *dist;
  NcmVector *x;
  gdouble rs_fiduc;
};

GType nc_data_bao_dmr_hr_get_type (void) G_GNUC_CONST;

NcDataBaoDMrHr *nc_data_bao_dmr_hr_new_from_file (const gchar *filename);
NcDataBaoDMrHr *nc_data_bao_dmr_hr_new_from_id (NcDistance *dist, NcDataBaoId id);
void nc_data_bao_dmr_hr_set_dist (NcDataBaoDMrHr *dmh, NcDistance *dist);

G_END_DECLS

#endif /* _NC_DATA_BAO_DMR_HR_H_ */

/***************************************************************************
 *            nc_data_bao_empirical_fit.h
 *
 *  Wed February 11 13:03:16 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_data_bao_empirical_fit.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_BAO_EMPIRICAL_FIT_H_
#define _NC_DATA_BAO_EMPIRICAL_FIT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_dist1d.h>
#include <numcosmo/math/ncm_stats_dist1d_spline.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/data/nc_data_bao.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_BAO_EMPIRICAL_FIT             (nc_data_bao_empirical_fit_get_type ())
#define NC_DATA_BAO_EMPIRICAL_FIT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_BAO_EMPIRICAL_FIT, NcDataBaoEmpiricalFit))
#define NC_DATA_BAO_EMPIRICAL_FIT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_BAO_EMPIRICAL_FIT, NcDataBaoEmpiricalFitClass))
#define NC_IS_DATA_BAO_EMPIRICAL_FIT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_BAO_EMPIRICAL_FIT))
#define NC_IS_DATA_BAO_EMPIRICAL_FIT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_BAO_EMPIRICAL_FIT))
#define NC_DATA_BAO_EMPIRICAL_FIT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_BAO_EMPIRICAL_FIT, NcDataBaoEmpiricalFitClass))

typedef struct _NcDataBaoEmpiricalFitClass NcDataBaoEmpiricalFitClass;
typedef struct _NcDataBaoEmpiricalFit NcDataBaoEmpiricalFit;

struct _NcDataBaoEmpiricalFitClass
{
  /*< private >*/
  NcmDataDist1dClass parent_class;
};

struct _NcDataBaoEmpiricalFit
{
  /*< private >*/
  NcmDataDist1d parent_instance;
  gdouble Dv_fiduc;
  gdouble rs_fiduc;
  gdouble z;
  NcmSpline *m2lnp;
  NcmStatsDist1d *p;
  gdouble p_mode;
  NcDistance *dist;
};

GType nc_data_bao_empirical_fit_get_type (void) G_GNUC_CONST;

NcDataBaoEmpiricalFit *nc_data_bao_empirical_fit_new (NcmSpline *m2lnp, gdouble Dv_fiduc, gdouble rs_fiduc, gdouble z);
NcDataBaoEmpiricalFit *nc_data_bao_empirical_fit_new_from_file (const gchar *filename);
NcDataBaoEmpiricalFit *nc_data_bao_empirical_fit_new_from_id (NcDistance *dist, NcDataBaoId id);

gdouble nc_data_bao_empirical_fit_get_mode (NcDataBaoEmpiricalFit *bao_ef);
gdouble nc_data_bao_empirical_fit_get_alpha (NcDataBaoEmpiricalFit *bao_ef, NcmMSet *mset);

void nc_data_bao_empirical_fit_set_dist (NcDataBaoEmpiricalFit *bao_ef, NcDistance *dist);

G_END_DECLS

#endif /* _NC_DATA_BAO_EMPIRICAL_FIT_H_ */

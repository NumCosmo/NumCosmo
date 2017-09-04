/***************************************************************************
 *            nc_data_bao_empirical_fit_2d.h
 *
 *  Fri September 1 14:39:23 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_bao_empirical_fit_2d.h
 * Copyright (C) 2017 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_DATA_BAO_EMPIRICAL_FIT_2D_H_
#define _NC_DATA_BAO_EMPIRICAL_FIT_2D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_dist2d.h>
#include <numcosmo/math/ncm_stats_dist2d_spline.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/data/nc_data_bao.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_BAO_EMPIRICAL_FIT_2D             (nc_data_bao_empirical_fit_2d_get_type ())
#define NC_DATA_BAO_EMPIRICAL_FIT_2D(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_BAO_EMPIRICAL_FIT_2D, NcDataBaoEmpiricalFit2d))
#define NC_DATA_BAO_EMPIRICAL_FIT_2D_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_BAO_EMPIRICAL_FIT_2D, NcDataBaoEmpiricalFit2dClass))
#define NC_IS_DATA_BAO_EMPIRICAL_FIT_2D(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_BAO_EMPIRICAL_FIT_2D))
#define NC_IS_DATA_BAO_EMPIRICAL_FIT_2D_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_BAO_EMPIRICAL_FIT_2D))
#define NC_DATA_BAO_EMPIRICAL_FIT_2D_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_BAO_EMPIRICAL_FIT_2D, NcDataBaoEmpiricalFit2dClass))

typedef struct _NcDataBaoEmpiricalFit2dClass NcDataBaoEmpiricalFit2dClass;
typedef struct _NcDataBaoEmpiricalFit2d NcDataBaoEmpiricalFit2d;

struct _NcDataBaoEmpiricalFit2dClass
{
  /*< private >*/
  NcmDataDist2dClass parent_class;
};

struct _NcDataBaoEmpiricalFit2d
{
  /*< private >*/
  NcmDataDist2d parent_instance;
  gdouble Dh_rd_fiduc;
  gdouble Dt_rd_fiduc;
  gdouble z;
  NcmSpline2d *m2lnp;
  NcmStatsDist2d *p;
  gdouble p_mode;
  NcDistance *dist;
};

GType nc_data_bao_empirical_fit_2d_get_type (void) G_GNUC_CONST;

NcDataBaoEmpiricalFit2d *nc_data_bao_empirical_fit_2d_new (NcmSpline2d *m2lnp, gdouble Dh_rd_fiduc, gdouble Dt_rd_fiduc, gdouble z);
NcDataBaoEmpiricalFit2d *nc_data_bao_empirical_fit_2d_new_from_file (const gchar *filename);
NcDataBaoEmpiricalFit2d *nc_data_bao_empirical_fit_2d_new_from_id (NcDistance *dist, NcDataBaoId id);

gdouble nc_data_bao_empirical_fit_2d_get_mode (NcDataBaoEmpiricalFit2d *bao_ef);
gdouble nc_data_bao_empirical_fit_2d_get_alpha_perpendicular (NcDataBaoEmpiricalFit2d *bao_ef, NcmMSet *mset);
gdouble nc_data_bao_empirical_fit_2d_get_alpha_parallel (NcDataBaoEmpiricalFit2d *bao_ef, NcmMSet *mset);

void nc_data_bao_empirical_fit_2d_set_dist (NcDataBaoEmpiricalFit2d *bao_ef, NcDistance *dist);

G_END_DECLS

#endif /* _NC_DATA_BAO_EMPIRICAL_FIT_2D_H_ */

/***************************************************************************
 *            nc_data_bao_empirical_fit.h
 *
 *  Wed February 11 13:03:16 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_bao_empirical_fit.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

G_DECLARE_FINAL_TYPE (NcDataBaoEmpiricalFit, nc_data_bao_empirical_fit, NC, DATA_BAO_EMPIRICAL_FIT, NcmDataDist1d);

NcDataBaoEmpiricalFit *nc_data_bao_empirical_fit_new (NcmSpline *m2lnp, gdouble Dv_fiduc, gdouble rs_fiduc, gdouble z);
NcDataBaoEmpiricalFit *nc_data_bao_empirical_fit_new_from_file (const gchar *filename);
NcDataBaoEmpiricalFit *nc_data_bao_empirical_fit_new_from_id (NcDistance *dist, NcDataBaoId id);

gdouble nc_data_bao_empirical_fit_get_mode (NcDataBaoEmpiricalFit *bao_ef);
gdouble nc_data_bao_empirical_fit_get_alpha (NcDataBaoEmpiricalFit *bao_ef, NcmMSet *mset);

void nc_data_bao_empirical_fit_set_dist (NcDataBaoEmpiricalFit *bao_ef, NcDistance *dist);

G_END_DECLS

#endif /* _NC_DATA_BAO_EMPIRICAL_FIT_H_ */


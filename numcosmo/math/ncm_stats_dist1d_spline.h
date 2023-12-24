/***************************************************************************
 *            ncm_stats_dist1d_spline.h
 *
 *  Thu February 12 16:50:52 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d_spline.h
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

#ifndef _NCM_STATS_DIST1D_SPLINE_H_
#define _NCM_STATS_DIST1D_SPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_stats_dist1d.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST1D_SPLINE (ncm_stats_dist1d_spline_get_type ())

G_DECLARE_FINAL_TYPE (NcmStatsDist1dSpline, ncm_stats_dist1d_spline, NCM, STATS_DIST1D_SPLINE, NcmStatsDist1d)


NcmStatsDist1dSpline *ncm_stats_dist1d_spline_new (NcmSpline * m2lnp);

G_END_DECLS

#endif /* _NCM_STATS_DIST1D_SPLINE_H_ */


/***************************************************************************
 *            ncm_stats_dist2d_spline.h
 *
 *  Sat July 22 22:31:17 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * ncm_stats_dist2d_spline.h
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

#ifndef _NCM_STATS_DIST2D_SPLINE_H_
#define _NCM_STATS_DIST2D_SPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline2d.h>
#include <numcosmo/math/ncm_stats_dist2d.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST2D_SPLINE (ncm_stats_dist2d_spline_get_type ())

G_DECLARE_FINAL_TYPE (NcmStatsDist2dSpline, ncm_stats_dist2d_spline, NCM, STATS_DIST2D_SPLINE, NcmStatsDist2d)

NcmStatsDist2dSpline *ncm_stats_dist2d_spline_new (NcmSpline2d * m2lnp);

G_END_DECLS

#endif /* _NCM_STATS_DIST2D_SPLINE_H_ */


/***************************************************************************
 *            ncm_stats_dist1d_spline.h
 *
 *  Thu February 12 16:50:52 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d_spline.h
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

#ifndef _NCM_STATS_DIST1D_SPLINE_H_
#define _NCM_STATS_DIST1D_SPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_stats_dist1d.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST1D_SPLINE             (ncm_stats_dist1d_spline_get_type ())
#define NCM_STATS_DIST1D_SPLINE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST1D_SPLINE, NcmStatsDist1dSpline))
#define NCM_STATS_DIST1D_SPLINE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST1D_SPLINE, NcmStatsDist1dSplineClass))
#define NCM_IS_STATS_DIST1D_SPLINE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST1D_SPLINE))
#define NCM_IS_STATS_DIST1D_SPLINE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST1D_SPLINE))
#define NCM_STATS_DIST1D_SPLINE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST1D_SPLINE, NcmStatsDist1dSplineClass))

typedef struct _NcmStatsDist1dSplineClass NcmStatsDist1dSplineClass;
typedef struct _NcmStatsDist1dSpline NcmStatsDist1dSpline;

struct _NcmStatsDist1dSplineClass
{
  /*< private >*/
  NcmStatsDist1dClass parent_class;
};

typedef struct _NcmStatsDist1dSplineTail
{
  gdouble xb;
  gdouble a;
  gdouble b;
  gdouble c;
  gdouble sigma;
} NcmStatsDist1dSplineTail;

struct _NcmStatsDist1dSpline
{
  /*< private >*/
  NcmStatsDist1d parent_instance;
  NcmSpline *m2lnp;
  gdouble tail_sigma;
  NcmStatsDist1dSplineTail left_tail;
  NcmStatsDist1dSplineTail right_tail;
};

GType ncm_stats_dist1d_spline_get_type (void) G_GNUC_CONST;

NcmStatsDist1dSpline *ncm_stats_dist1d_spline_new (NcmSpline *m2lnp);

G_END_DECLS

#endif /* _NCM_STATS_DIST1D_SPLINE_H_ */


/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_stats_dist_nd_kde_studentt.h
 *
 *  Wed November 07 17:41:38 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd_kde_studentt.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_STATS_DIST_ND_KDE_STUDENTT_H_
#define _NCM_STATS_DIST_ND_KDE_STUDENTT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_stats_dist_nd.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_ND_KDE_STUDENTT             (ncm_stats_dist_nd_kde_studentt_get_type ())
#define NCM_STATS_DIST_ND_KDE_STUDENTT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST_ND_KDE_STUDENTT, NcmStatsDistNdKDEStudentt))
#define NCM_STATS_DIST_ND_KDE_STUDENTT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST_ND_KDE_STUDENTT, NcmStatsDistNdKDEStudenttClass))
#define NCM_IS_STATS_DIST_ND_KDE_STUDENTT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST_ND_KDE_STUDENTT))
#define NCM_IS_STATS_DIST_ND_KDE_STUDENTT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST_ND_KDE_STUDENTT))
#define NCM_STATS_DIST_ND_KDE_STUDENTT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST_ND_KDE_STUDENTT, NcmStatsDistNdKDEStudenttClass))

typedef struct _NcmStatsDistNdKDEStudenttClass NcmStatsDistNdKDEStudenttClass;
typedef struct _NcmStatsDistNdKDEStudentt NcmStatsDistNdKDEStudentt;
typedef struct _NcmStatsDistNdKDEStudenttPrivate NcmStatsDistNdKDEStudenttPrivate;


struct _NcmStatsDistNdKDEStudenttClass
{
  NcmStatsDistNdClass parent_class;
};

struct _NcmStatsDistNdKDEStudentt
{
  NcmStatsDistNd parent_instance;

  NcmStatsDistNdKDEStudenttPrivate *priv;
};

GType ncm_stats_dist_nd_kde_studentt_get_type (void) G_GNUC_CONST;

NcmStatsDistNdKDEStudentt *ncm_stats_dist_nd_kde_studentt_new (const guint dim, const gboolean LOOCV, const guint freedom);

NcmStatsDistNdKDEStudentt *ncm_stats_dist_nd_kde_studentt_ref (NcmStatsDistNdKDEStudentt *dndg);
void ncm_stats_dist_nd_kde_studentt_free (NcmStatsDistNdKDEStudentt *dndg);
void ncm_stats_dist_nd_kde_studentt_clear (NcmStatsDistNdKDEStudentt **dndg);

void ncm_stats_dist_nd_kde_studentt_set_over_smooth (NcmStatsDistNdKDEStudentt *dndg, const gdouble over_smooth);
gdouble ncm_stats_dist_nd_kde_studentt_get_over_smooth (NcmStatsDistNdKDEStudentt *dndg);

void ncm_stats_dist_nd_kde_studentt_set_LOOCV_bandwidth_adj (NcmStatsDistNdKDEStudentt *dndg, gboolean LOOCV);
gboolean ncm_stats_dist_nd_kde_studentt_get_LOOCV_bandwidth_adj (NcmStatsDistNdKDEStudentt *dndg);

void ncm_stats_dist_nd_kde_studentt_add_obs_weight (NcmStatsDistNdKDEStudentt *dndg, NcmVector *y, const gdouble w);
void ncm_stats_dist_nd_kde_studentt_add_obs (NcmStatsDistNdKDEStudentt *dndg, NcmVector *y);

G_END_DECLS

#endif /* _NCM_STATS_DIST_ND_KDE_STUDENTT_H_ */

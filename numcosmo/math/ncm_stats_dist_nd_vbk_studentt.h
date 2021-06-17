/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_nd_vbk_studentt.h
 *
 *  Wed November 07 17:41:38 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd_vbk_studentt.h
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

#ifndef _NCM_STATS_DIST_ND_VBK_STUDENTT_H_
#define _NCM_STATS_DIST_ND_VBK_STUDENTT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_stats_dist_nd_vbk.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_ND_VBK_STUDENTT             (ncm_stats_dist_nd_vbk_studentt_get_type ())
#define NCM_STATS_DIST_ND_VBK_STUDENTT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST_ND_VBK_STUDENTT, NcmStatsDistNdVBKStudentt))
#define NCM_STATS_DIST_ND_VBK_STUDENTT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST_ND_VBK_STUDENTT, NcmStatsDistNdVBKStudenttClass))
#define NCM_IS_STATS_DIST_ND_VBK_STUDENTT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST_ND_VBK_STUDENTT))
#define NCM_IS_STATS_DIST_ND_VBK_STUDENTT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST_ND_VBK_STUDENTT))
#define NCM_STATS_DIST_ND_VBK_STUDENTT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST_ND_VBK_STUDENTT, NcmStatsDistNdVBKStudenttClass))

typedef struct _NcmStatsDistNdVBKStudenttClass NcmStatsDistNdVBKStudenttClass;
typedef struct _NcmStatsDistNdVBKStudentt NcmStatsDistNdVBKStudentt;
typedef struct _NcmStatsDistNdVBKStudenttPrivate NcmStatsDistNdVBKStudenttPrivate;


struct _NcmStatsDistNdVBKStudenttClass
{
  NcmStatsDistNdVBKClass parent_class;
};

struct _NcmStatsDistNdVBKStudentt
{
  NcmStatsDistNdVBK parent_instance;
  
  NcmStatsDistNdVBKStudenttPrivate *priv;
};

GType ncm_stats_dist_nd_vbk_studentt_get_type (void) G_GNUC_CONST;

NcmStatsDistNdVBKStudentt *ncm_stats_dist_nd_vbk_studentt_new (const guint dim, const NcmStatsDistNdVBKCV cv_type, const gdouble nu);

NcmStatsDistNdVBKStudentt *ncm_stats_dist_nd_vbk_studentt_ref (NcmStatsDistNdVBKStudentt *dndt);
void ncm_stats_dist_nd_vbk_studentt_free (NcmStatsDistNdVBKStudentt *dndt);
void ncm_stats_dist_nd_vbk_studentt_clear (NcmStatsDistNdVBKStudentt **dndt);

void ncm_stats_dist_nd_vbk_studentt_set_nu (NcmStatsDistNdVBKStudentt *dndt, const gdouble nu);
gdouble ncm_stats_dist_nd_vbk_studentt_get_nu (NcmStatsDistNdVBKStudentt *dndt);

G_END_DECLS

#endif /* _NCM_STATS_DIST_ND_VBK_STUDENTT_H_ */


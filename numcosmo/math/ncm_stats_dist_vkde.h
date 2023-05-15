/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_vkde.h
 *
 *  Wed November 07 16:02:25 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_vkde.h
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_STATS_DIST_VKDE_H_
#define _NCM_STATS_DIST_VKDE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_stats_vec.h>
#include <numcosmo/math/ncm_stats_dist_kde.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_VKDE             (ncm_stats_dist_vkde_get_type ())
#define NCM_STATS_DIST_VKDE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST_VKDE, NcmStatsDistVKDE))
#define NCM_STATS_DIST_VKDE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST_VKDE, NcmStatsDistVKDEClass))
#define NCM_IS_STATS_DIST_VKDE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST_VKDE))
#define NCM_IS_STATS_DIST_VKDE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST_VKDE))
#define NCM_STATS_DIST_VKDE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST_VKDE, NcmStatsDistVKDEClass))

typedef struct _NcmStatsDistVKDEClass NcmStatsDistVKDEClass;
typedef struct _NcmStatsDistVKDE NcmStatsDistVKDE;
typedef struct _NcmStatsDistVKDEPrivate NcmStatsDistVKDEPrivate;

struct _NcmStatsDistVKDEClass
{
  /*< private >*/
  NcmStatsDistKDEClass parent_class;
};

struct _NcmStatsDistVKDE
{
  /*< private >*/
  NcmStatsDistKDE parent_instance;
  NcmStatsDistVKDEPrivate *priv;
};

GType ncm_stats_dist_vkde_get_type (void) G_GNUC_CONST;

NcmStatsDistVKDE *ncm_stats_dist_vkde_new (NcmStatsDistKernel *sdk, NcmStatsDistCV CV_type);
NcmStatsDistVKDE *ncm_stats_dist_vkde_ref (NcmStatsDistVKDE *sdvkde);
void ncm_stats_dist_vkde_free (NcmStatsDistVKDE *sdvkde);
void ncm_stats_dist_vkde_clear (NcmStatsDistVKDE **sdvkde);

void ncm_stats_dist_vkde_set_local_frac (NcmStatsDistVKDE *sdvkde, const gdouble local_frac);
gdouble ncm_stats_dist_vkde_get_local_frac (NcmStatsDistVKDE *sdvkde);

void ncm_stats_dist_vkde_set_use_rot_href (NcmStatsDistVKDE *sdvkde, const gboolean use_rot_href);
gboolean ncm_stats_dist_vkde_get_use_rot_href (NcmStatsDistVKDE *sdvkde);

void ncm_stats_dist_vkde_set_use_threads (NcmStatsDistVKDE *sdvkde, const gboolean use_threads);
gboolean ncm_stats_dist_vkde_get_use_threads (NcmStatsDistVKDE *sdvkde);

G_END_DECLS

#endif /* _NCM_STATS_DIST_ND_H_ */


/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_stats_dist_nd.h
 *
 *  Wed November 07 16:02:25 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd.h
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

#ifndef _NCM_STATS_DIST_ND_H_
#define _NCM_STATS_DIST_ND_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_ND             (ncm_stats_dist_nd_get_type ())
#define NCM_STATS_DIST_ND(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST_ND, NcmStatsDistNd))
#define NCM_STATS_DIST_ND_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST_ND, NcmStatsDistNdClass))
#define NCM_IS_STATS_DIST_ND(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST_ND))
#define NCM_IS_STATS_DIST_ND_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST_ND))
#define NCM_STATS_DIST_ND_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST_ND, NcmStatsDistNdClass))

typedef struct _NcmStatsDistNdClass NcmStatsDistNdClass;
typedef struct _NcmStatsDistNd NcmStatsDistNd;
typedef struct _NcmStatsDistNdPrivate NcmStatsDistNdPrivate;

struct _NcmStatsDistNdClass
{
  /*< private >*/ 
  GObjectClass parent_class;
  void (*prepare) (NcmStatsDistNd *dnd);
  void (*prepare_interp) (NcmStatsDistNd *dnd, NcmVector *m2lnp, gboolean normalized);
  void (*set_dim) (NcmStatsDistNd *dnd, const guint dim);
  gdouble (*eval) (NcmStatsDistNd *dnd, NcmVector *x);
  gdouble (*eval_m2lnp) (NcmStatsDistNd *dnd, NcmVector *x);
  void (*sample) (NcmStatsDistNd *dnd, NcmVector *x, NcmRNG *rng);
};

struct _NcmStatsDistNd
{
  /*< private >*/
  GObject parent_instance;
  NcmStatsDistNdPrivate *priv;
};

GType ncm_stats_dist_nd_get_type (void) G_GNUC_CONST;

NcmStatsDistNd *ncm_stats_dist_nd_ref (NcmStatsDistNd *dnd);
void ncm_stats_dist_nd_free (NcmStatsDistNd *dnd);
void ncm_stats_dist_nd_clear (NcmStatsDistNd **dnd);

guint ncm_stats_dist_nd_get_dim (NcmStatsDistNd *dnd);

void ncm_stats_dist_nd_prepare (NcmStatsDistNd *dnd);
void ncm_stats_dist_nd_prepare_interp (NcmStatsDistNd *dnd, NcmVector *m2lnp, gboolean normalized);
gdouble ncm_stats_dist_nd_eval (NcmStatsDistNd *dnd, NcmVector *x);
gdouble ncm_stats_dist_nd_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *x);
void ncm_stats_dist_nd_sample (NcmStatsDistNd *dnd, NcmVector *x, NcmRNG *rng);

G_END_DECLS

#endif /* _NCM_STATS_DIST_ND_H_ */
